FROM maven:3-jdk-11 as builder

COPY . /usr/src/app
WORKDIR /usr/src/app

RUN mvn clean compile assembly:single

FROM openjdk:11-jre-slim
WORKDIR /usr/src/app

COPY ./template/template_v9.html .

COPY --from=builder /usr/src/app/target/multi-tau-correlator-acf-1.0.0-jar-with-dependencies.jar .
CMD [ "java", "-jar", "multi-tau-correlator-acf-1.0.0-jar-with-dependencies.jar" ]
